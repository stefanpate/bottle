# Biological computer aided synthesis planning

## Overview

This is a computational pipeline designed for biological synthesis planning. The system performs network expansions using Pickaxe, discovers synthetic pathways, processes results, and visualizes results in an interactive app.

### Key Features

- 🧪 **Network Expansion**: Automated chemical reaction network generation using Pickaxe
- 🔍 **Pathway Discovery**: Advanced algorithms for finding optimal synthetic routes
- 📊 **Structural Analysis**: Machine learning-based feasibility scoring for predicted reactions
- ⚡ **Thermodynamic Calculations**: Integration with eQuilibrator for thermodynamic analysis
- 📈 **Visualization**: Interactive pathway viewer and reaction drawing tools

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

Additionally, create a `.env` file at the project root by copying `.ENV_TEMPLATE`:

```bash
cp .ENV_TEMPLATE .env
```

Then edit `.env` and set both variables to absolute paths on your machine:

- `CASP_STUDY_ROOT` — directory containing your processed CASP study subdirectories (used by the path viewer and baked into the Docker image at build time).
- `KNOWN_BIOCHEM_ROOT` — directory containing the known biochemistry tables (typically `artifacts/known` in this repo).

These are loaded automatically by `poe` tasks (via `envfile = ".env"` in `pyproject.toml`).

#### 2. Download Equilibrator Database

The below command will download the eQuilibrator database to a canonical cache location depending on your operating system, e.g., linux: `~/.cache/equilibrator`.

```python
import equilibrator_assets.local_compound_cache
```

After the download is complete, optionally move the downloaded subdir and configure the     `equilibrator_cache` path in your `filepaths.yaml` appropriately.

## 📁 Project Structure

1. **artifacts** - Ancillary small data files
2. **conf** - Configuration files for the pipeline
3. **logs** - Log files and error tracking
4. **scripts** - Main executable scripts
5. **notebooks** - Jupyter notebooks for analysis
6. **src** - Source packages and modules
7. **path_viewer** - Visualization tools
8. **tests** - Test suite

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
python linear_pathfinding.py forward_expansion=my/fwd/subdir retro_expansion=my/retro/subdir casp_study=my_casp_study
```

> **Note:** `retrosynthesis.py` utilizes a network of known enzymatic reactions in addition to predicted reactions from the expansion to find synthetic paths.

**2. Reaction structural analysis**

Looks up known reactions similar to the predicted reactions and collects their enzymes. Feasibility of the predicted reactions is scored with a [machine learning model](https://github.com/tyo-nu/DORA_XGB).

**Example usage:**
```bash
python analyze_structures.py casp_study=my_casp_study
```

**3. Thermodynamic calculations (Optional)**

Calculates thermodynamic values (e.g., ΔG) for predicted pathways using [eQuilibrator](https://equilibrator.readthedocs.io/en/latest/local_cache.html).

**Example usage:**
```bash
python analyze_thermo.py casp_study=my_casp_study
```

**4. Reaction drawing**

Draws reactions to SVG format.

**Example usage:**
```bash
python draw_reactions.py casp_study=my_casp_study
```

### 🖥️ Visualizing Processed Synthesis Paths

Processed paths can be visualized using the interactive path viewer app. To launch it locally:

```bash
streamlit run path_viewer/list_view.py -- --casp-study my_casp_study
```

Where `my_casp_study` is the name of the subdirectory in your processed data location that contains the processed parquet and SVG files.

### 🐳 Containerizing the Path Viewer

The path viewer can be built and run as a Docker image using [poe the poet](https://poethepoet.natn.io/) tasks defined in `pyproject.toml`. `poe` lives in the `dev` dependency group, so install it with:

```bash
poetry install --with dev
```

Make sure your `.env` file is set up (see [Configuration Setup](#️-configuration-setup)) — `CASP_STUDY_ROOT` and `KNOWN_BIOCHEM_ROOT` are read at build time and baked into the image.

The main tasks:

```bash
poetry run poe build-image   # Build the image, baking in CASP and known data
poetry run poe run-docker    # Run the container, exposing the viewer on :8501
poetry run poe push-image    # Push the image to the registry
```

Once running, the viewer is available at <http://localhost:8501>.

## 🔧 Using Other Network Expansion Software

Users may use different network expansion software upstream of the BOTTLE post-processing pipeline, provided they output data adhering to the interim data schemas defined in `src/schemas.py`.

Each schema specifies required columns that must be present in the output files. The expected columns for each schema are:

### 📋 Data Schema Requirements

#### 1. Expansion Reactions Schema

| Column Name       | Description                                                |
|-------------------|------------------------------------------------------------|
| `am_smarts`       | Atom-mapped reaction SMARTS                                |
| `rule_name`       | Name of rule that generated the reaction                  |
| `rule_set`        | Name of rule set                                            |
| `rule_template`   | SMARTS-encoded rule template                               |

#### 2. Compounds Schema

| Column Name | Description                                                     |
|-------------|-----------------------------------------------------------------|
| `id`        | Unique identifier for the compound                              |
| `smiles`    | SMILES string representation of the compound                   |
| `type`      | Compound type: `'source'`, `'target'`, `'known'`, `'helper'`, or `'checkpoint'` |
| `name`      | Compound name (optional)                                       |

---
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**BOTTLE** (Biological Computer Aided Synthesis Planning) is a computational pipeline for discovering biological synthesis routes. It expands chemical reaction networks, finds synthetic pathways, scores reaction feasibility (structurally and thermodynamically), and visualizes results in a Streamlit UI.

## Commands

### Installation
```bash
poetry install           # Install all dependencies
poetry install --with dev  # Include dev tools (pytest, memory-profiler, findimports)
poetry install --with lab  # Include Jupyter lab extras
```

### Tests
```bash
poetry run pytest --cov=src          # Run all tests with coverage
poetry run pytest tests/test_smarts.py  # Run a single test file
poetry run poe test                  # Alias via poe
```

### Pipeline Execution (run in order)
```bash
# 1. Expand chemical network
python scripts/expand.py rules=mechinformed_rules starters=my_starters targets=my_targets

# 2a. Find linear paths
python scripts/linear_pathfinding.py casp_study=my_study expansion=path/to/expansion

# 2b. OR: Retrosynthetic path finding
python scripts/retrosynthesis.py casp_study=my_study

# 3. Score reaction feasibility
python scripts/analyze_structures.py casp_study=my_study

# 4. Thermodynamic analysis (optional)
python scripts/analyze_thermo.py casp_study=my_study

# 5. Generate reaction SVGs
python scripts/draw_reactions.py casp_study=my_study

# 6. Launch interactive viewer
streamlit run path_viewer/list_view.py -- --casp-study my_study
```

### Docker
```bash
poetry run poe build-image    # Build amd64 Docker image
poetry run poe run-docker     # Run in Docker
```

## Architecture

### Pipeline Data Flow

The pipeline stages communicate via parquet files. Each stage reads upstream parquets and writes its own:

```
expand.py
  → raw_data/expansions/{expansion}/
      reactions.parquet     (am_smarts, rule_name, rule_set)
      compounds.parquet     (id, smiles, type, name)

linear_pathfinding.py / retrosynthesis.py
  → interim_data/{casp_study}/
      predicted_reactions.parquet
      paths.parquet              (path_id, rxn_id, generation, rxn_type)
      path_stats.parquet         (id, starters, targets, mdf, feasibility_frac, ...)
      compounds.parquet

analyze_structures.py
  → predicted_reactions.parquet  [adds: dxgb_label, rxn_sims, analogue_ids]

analyze_thermo.py
  → path_stats.parquet           [adds: mdf, dg_opt, dg_err]

draw_reactions.py
  → processed_data/{casp_study}/svgs/

path_viewer/list_view.py        (reads all processed outputs)
```

### Configuration

All scripts use **Hydra** for configuration management. YAML configs live in `conf/`. File paths are configured separately via `conf/filepaths/filepaths.yaml` (created from `FILEPATHS_TEMPLATE.yaml`). Override any config value on the command line with `key=value`.

### Core Modules (`src/`)

**`network.py`** — Central data structure. `ReactionNetwork` extends `networkx.MultiDiGraph` where nodes are compounds and edges are reactions. Key concepts:
- `pnmc` (product normalized mass contribution) and `rnmc` (reaction normalized mass contribution) quantify atom flow through edges
- `prune()` removes low-contribution edges by PNMC/RNMC thresholds
- `enumerate_synthetic_trees()` generates feasible pathways as `SyntheticTree` objects via depth-first enumeration
- `add_reaction()` expects atom-mapped SMARTS strings

**`schemas.py`** — Polars schemas for all parquet files. Reference here when reading/writing pipeline data.

**`post_processing.py`** — `PathWrangler` class loads and filters processed outputs. `pick_constraints_for_MDF()` generates eQuilibrator concentration constraints per cofactor system (ATP/ADP, NAD(P)H, etc.).

**`chem_draw.py`** — RDKit + svgutils rendering of reactions and molecules to SVG.

### Key Algorithms

**Network expansion** (`expand.py`): Uses Pickaxe (minedatabase) to apply reaction rules to starter compounds, optionally sampling by Tanimoto similarity to targets.

**Linear pathfinding** (`linear_pathfinding.py`): Constructs a `ReactionNetwork`, prunes by mass contribution thresholds, then BFS/DFS finds half-paths (forward from starters, retro from targets) and combines them.

**Retrosynthesis** (`retrosynthesis.py`): Merges predicted reactions into a known enzymatic reaction network (JSON), then uses `enumerate_synthetic_trees()` to recursively build all possible retrosynthetic trees.

**Structure analysis** (`analyze_structures.py`): Computes Reaction Center Morgan Fingerprints (RCMFP) for predicted reactions, calculates Tanimoto similarity to a known reaction database, and scores feasibility with DORA-XGB.

**Thermodynamic analysis** (`analyze_thermo.py`): Uses eQuilibrator to calculate ΔG values and formulates max-min driving force (MDF) as a cvxpy optimization with cofactor-specific concentration constraints.

### Compound Types

Compounds in the network have a `type` field:
- `source` — starting materials (no synthesis required)
- `target` — desired products
- `helper` — cofactors/coreactants
- `known` — in the known reaction database
- `intermediate` — network-generated intermediates
- `checkpoint` — user-defined intermediate waypoints

### External Requirements

- **ChemAxon `cxcalc`**: Required for pKa calculations in `analyze_thermo.py`
- **eQuilibrator database**: Must be configured for thermodynamic analysis
- **DORA-XGB model**: Required for feasibility scoring in `analyze_structures.py`

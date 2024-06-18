# Computationally-assisted synthesis planning for recyclable-by-design polymers

## Description
Runs network expansions using Pickaxe, finds paths, post-processes, and renders spreadsheet + pdf deliverables.

## Installation. 

## Directories
1. artifacts (stores results)
2. data (input data: known reactions, mappings, starters, targets, rules)
3. logs
4. scripts
5. notebooks
6. src (packages and modules)

## Usage
- run_pickaxe.py does a network expansion from starters to targets. 
- prune_pickaxe.py takes a pickaxe output, finds paths from starters to targets, saves a pruned pickaxe object and a ProcessedExpansion object
- [Placeholder for thermodynamics. Currently calculations done outside this directory]
- analyze_paths.py does reaction-center-maximum-common-substructure calculations and enters thermo info into ProcessedExpansion object
- vis_pathways.ipynb filters & sorts paths in ProcessedExpansion object and renders pdfs of reaction formula and spreadsheets w/ further info

## Dev Notes
1. 

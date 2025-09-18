#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10GB
#SBATCH -t 1:00:00
#SBATCH --job-name="lin_path"
#SBATCH --output=/home/spn1560/bottle/logs/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
script=/home/spn1560/bottle/scripts/linear_pathfinding.py
expansion=3_steps_bottle_targets_25_to_None_rules_mechinferred_dt_04_rules_w_coreactants_co_metacyc_coreactants_sampled_False_pruned_False_aplusb_False
casp_study=perf
max_depth=3

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script max_depth=$max_depth expansion=$expansion casp_study=$casp_study
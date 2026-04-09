#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=32GB
#SBATCH -t 4:00:00
#SBATCH --job-name="retro"
#SBATCH --output=/home/spn1560/bottle/logs/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
script=/home/spn1560/bottle/scripts/retrosynthesis.py
expansion=5_steps_ccm_aa_to_bottle_targets_25_combo_mechinferred_dt_13_and_mechinferred_dt_04_rules
casp_study=branched_test
max_depth=5
max_leaves=3
processes=10

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script max_depth=$max_depth max_leaves=$max_leaves expansion=$expansion casp_study=$casp_study processes=$processes